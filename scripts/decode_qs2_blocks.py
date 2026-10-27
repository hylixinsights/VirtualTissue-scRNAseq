"""Decode qs2 format v1 to an R serialization stream with bounded buffers.

Format references: qs2 0.3.1's file_headers.h, io_common.h, zstd_module.h and
xxhash_module.h. This independent decoder writes one <=1 MiB block at a time.
It validates the internal XXH3 checksum before promoting output. It does not
deserialize R objects. Validate byte equality with qs2::qs_to_rds fixtures.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'.pydeps-expansion'))
import numpy as np
import xxhash
import zstandard


def decode(source, target):
    if target.exists():
        raise FileExistsError(target)
    temporary = target.with_suffix(target.suffix+'.partial')
    total = blocks = shuffled_blocks = 0
    check = xxhash.xxh3_64()
    output_hash = hashlib.sha256()
    decompressor = zstandard.ZstdDecompressor()
    with source.open('rb') as inp, temporary.open('wb') as out:
        header = inp.read(24)
        if len(header) != 24 or header[:4] != bytes([11,14,10,193]):
            raise ValueError('Not a qs2 file')
        if header[4] != 1 or header[5] != 1 or header[6] not in [1,2] or header[7] not in [0,1]:
            raise ValueError('Unsupported qs2 version/compression/endian/shuffle')
        endian = 'little' if header[6] == 2 else 'big'
        stored = int.from_bytes(header[16:24], endian)
        if not stored:
            raise ValueError('Missing source internal checksum')
        while prefix := inp.read(4):
            if len(prefix) != 4:
                raise ValueError('Truncated block-size prefix')
            size_flags = int.from_bytes(prefix, endian)
            shuffled = bool(size_flags & (1 << 31))
            size = size_flags & 0x7fffffff
            if not 0 < size <= 1_100_000 or (shuffled and not header[7]):
                raise ValueError('Invalid block size or shuffle flag')
            compressed = inp.read(size)
            if len(compressed) != size:
                raise ValueError('Truncated compressed block')
            check.update(prefix); check.update(compressed)
            content_size = zstandard.frame_content_size(compressed)
            if content_size not in [zstandard.CONTENTSIZE_UNKNOWN] and not 0 <= content_size <= 1_048_576:
                raise ValueError('Uncompressed block exceeds 1 MiB')
            block = decompressor.decompress(compressed, max_output_size=1_048_576)
            if len(block) > 1_048_576:
                raise ValueError('Uncompressed block exceeds 1 MiB')
            if shuffled:
                length = len(block) - len(block) % 8
                raw = np.frombuffer(block[:length], dtype=np.uint8)
                block = raw.reshape(8, length//8).T.tobytes() + block[length:]
                shuffled_blocks += 1
            out.write(block); output_hash.update(block)
            blocks += 1; total += len(block)
            if blocks % 1000 == 0:
                print(json.dumps({'blocks':blocks,'decoded_bytes':total}),flush=True)
        actual = check.intdigest() or 1
        if actual != stored:
            raise ValueError('Internal source XXH3 checksum mismatch')
    temporary.rename(target)
    return {'blocks':blocks,'shuffled_blocks':shuffled_blocks,'decoded_bytes':total,
            'internal_xxh3':str(actual),'internal_checksum_verified':True,
            'output_sha256':output_hash.hexdigest()}


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('input',type=Path); cli.add_argument('output',type=Path)
    args = cli.parse_args()
    audit = decode(args.input,args.output)
    args.output.with_suffix('.integrity.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit),flush=True)
