"""Compare the independent bounded decoder with the official R implementation."""
import gzip
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

from decode_qs2_blocks import ROOT, decode, zstandard

with tempfile.TemporaryDirectory(prefix='vtsc-qs2-blocks-') as temporary:
    directory = Path(temporary)
    code = r'''
    .libPaths(c('.Rlib-expansion','.Rlib',.Library))
    d <- commandArgs(trailingOnly=TRUE)[[1]]
    set.seed(101)
    x <- list(matrix=matrix(rep(c(1.1,2.2,3.3,4.4),1000000),ncol=100),
              random=runif(20000),tail='unaligned_tail_bytes',
              metadata=data.frame(patient=c('a','b'),lineage=factor(c('CD8','CD4'))))
    for (shuffle in c(FALSE,TRUE)) {
      name <- as.character(shuffle)
      qs2::qs_save(x,file.path(d,paste0(name,'.qs2')),shuffle=shuffle)
      qs2::qs_to_rds(file.path(d,paste0(name,'.qs2')),file.path(d,paste0(name,'.rds')))
    }
    '''
    script = directory/'fixture.R';script.write_text(code)
    subprocess.run(['Rscript',str(script),str(directory)],cwd=ROOT,check=True)
    checks = {}
    for flag in ['FALSE','TRUE']:
        source = directory/f'{flag}.qs2';target=directory/f'{flag}.decoded'
        audit = decode(source,target)
        with gzip.open(directory/f'{flag}.rds','rb') as handle:
            expected = hashlib.sha256(handle.read()).hexdigest()
        assert audit['output_sha256'] == expected
        if flag == 'TRUE':
            assert audit['shuffled_blocks'] > 0
        checks[flag] = {**audit,'byte_equivalent_to_official_R':True}
    corrupt = bytearray((directory/'FALSE.qs2').read_bytes())
    corrupt[16] ^= 1  # corrupt the stored hash, preserving all compressed frames
    (directory/'bad.qs2').write_bytes(corrupt)
    try:
        decode(directory/'bad.qs2',directory/'bad.decoded')
    except ValueError as e:
        assert 'checksum mismatch' in str(e)
    else:
        raise AssertionError('Corrupt source checksum accepted')
    assert not (directory/'bad.decoded').exists()
    (ROOT/'reports/evidence/cohort-expansion/block-decoder-validation.json').write_text(
        json.dumps({'synthetic_only':True,'checks':checks,'corruption_rejected':True},indent=2)+'\n')
    print('PASS: exact official-R bytes with/without shuffling; corrupt checksum rejected.')
