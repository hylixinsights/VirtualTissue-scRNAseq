"""Fallback: parse an RDS stream while discarding oversized numeric arrays.

This is ONLY a metadata extractor, never a matrix exporter. Dropped arrays are
explicitly logged; the full Seurat object must never be converted or saved.
The metadata subtree is rejected if it contains any discarded numeric array.
Requires the separately pinned rdata extension. Test on a real synthetic Seurat
serialization before using this fallback on source data.
"""
import argparse
import gzip
import hashlib
import json
import resource
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / '.pydeps-expansion'))
import numpy as np
import rdata
from rdata.parser._parser import DEFAULT_ALTREP_MAP, Parser, RData, RObjectType
from rdata.parser._binary import ParserBinary


class StreamingParser(ParserBinary):
    def __init__(self, handle, marker, max_elements):
        Parser.__init__(self, expand_altrep=True, altrep_constructor_dict=DEFAULT_ALTREP_MAP)
        self.file = handle
        self.byteorder = '>' if marker == b'X\n' else self._detect_byteorder(memoryview(handle.peek(4)[:4]))
        self.max_elements = max_elements
        self.discarded = []
        self.discarded_values = set()
        self.discarded_array_refs = []  # prevent object-id reuse after dtype conversion
        self.parsed_nodes = 0

    def _parse_array_values(self, dtype, length):
        dtype = np.dtype(dtype)
        if length > self.max_elements:
            start = self.file.tell()
            self.file.seek(length * dtype.itemsize, 1)
            assert self.file.tell() == start + length * dtype.itemsize
            array = np.empty(0, dtype=dtype)
            self.discarded_values.add(id(array))
            self.discarded_array_refs.append(array)
            self.discarded.append({'uncompressed_offset': start, 'elements': length, 'dtype': str(dtype)})
            return array
        return super()._parse_array_values(dtype, length)

    def parse_R_object(self, *args, **kwargs):
        self.parsed_nodes += 1
        if self.parsed_nodes % 100_000 == 0:
            peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            peak_bytes = peak if sys.platform == 'darwin' else peak * 1024
            if peak_bytes > 7 * 1024**3:
                raise MemoryError('Metadata parser exceeded its 7 GiB process budget')
            print(json.dumps({'parsed_nodes': self.parsed_nodes,
                              'uncompressed_bytes_read': self.file.tell(),
                              'peak_rss_bytes': peak_bytes}), flush=True)
        before = len(self.discarded)
        node = super().parse_R_object(*args, **kwargs)
        if len(self.discarded) > before and node.info.type in {
            RObjectType.LGL, RObjectType.INT, RObjectType.REAL, RObjectType.CPLX,
        }:
            # Logical conversion or NA masking may return a different array.
            self.discarded_values.add(id(node.value))
            self.discarded_array_refs.append(node.value)
        return node


def dereference(node):
    while node is not None and node.info.type == RObjectType.REF:
        node = node.referenced_object
    return node


def attributes(node):
    out = {}
    pair = node.attributes
    while pair is not None and pair.info.type not in {RObjectType.NIL, RObjectType.NILVALUE}:
        assert pair.info.type == RObjectType.LIST
        tag = dereference(pair.tag)
        assert tag.info.type == RObjectType.SYM
        key = dereference(tag.value).value.decode('utf-8')
        value, pair = pair.value
        out[key] = value
    return out


def contains_discarded(node, discarded_ids, visited=None):
    if visited is None:
        visited = set()
    if node is None or id(node) in visited:
        return False
    visited.add(id(node))
    if id(node.value) in discarded_ids:
        return True
    if isinstance(node.value, (tuple, list)):
        if any(contains_discarded(x, discarded_ids, visited) for x in node.value):
            return True
    return any(contains_discarded(x, discarded_ids, visited)
               for x in [node.tag, node.attributes, node.referenced_object])


def extract(path, max_elements=2_000_000):
    import io
    with path.open('rb') as source:
        gzipped = source.read(2) == b'\x1f\x8b'
    opener = gzip.open if gzipped else open
    with opener(path, 'rb') as source:
        handle = io.BufferedReader(source)
        marker = handle.read(2)
        assert marker in {b'B\n', b'X\n'}, marker
        parser = StreamingParser(handle, marker, max_elements)
        parsed = parser.parse_all()
        parser.check_complete()
    assert parsed.object.info.type == RObjectType.S4
    slots = attributes(parsed.object)
    metadata = slots['meta.data']
    assert not contains_discarded(metadata, parser.discarded_values), 'Metadata contains a discarded array'
    result = rdata.conversion.convert(RData(parsed.versions, parsed.extra, metadata))
    assert result.index.is_unique
    assert len(result) > 0
    return result, {'rdata_version': rdata.__version__, 'slot_names': list(slots),
                    'discarded_numeric_arrays': parser.discarded,
                    'metadata_rows': len(result), 'metadata_columns': list(result.columns),
                    'metadata_numeric_arrays_discarded': False,
                    'matrices_exported': False, 'state_changes_computed': False}


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('input', type=Path)
    cli.add_argument('output', type=Path)
    cli.add_argument('--max-elements', type=int, default=2_000_000)
    args = cli.parse_args()
    result, audit = extract(args.input, args.max_elements)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(args.output)
    args.output.with_suffix('.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps({'metadata_rows': len(result), 'metadata_columns': list(result.columns),
                      'discarded_numeric_arrays': len(audit['discarded_numeric_arrays'])}))
