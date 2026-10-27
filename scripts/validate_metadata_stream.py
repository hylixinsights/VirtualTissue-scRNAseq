"""Independent R-to-Python metadata comparison and rejection test on fixtures."""
import json
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from extract_seurat_metadata_stream import ROOT, extract

with tempfile.TemporaryDirectory(prefix='vtsc-stream-') as temporary:
    directory = Path(temporary)
    r_source = r'''
    .libPaths(c('.Rlib-expansion','.Rlib',.Library))
    suppressPackageStartupMessages(library(SeuratObject))
    suppressPackageStartupMessages(library(Matrix))
    args <- commandArgs(trailingOnly=TRUE); d <- args[[1]]
    x <- Matrix(matrix(rep(c(1L,2L,3L),40000),nrow=300),sparse=TRUE)
    rownames(x) <- paste0('gene',seq_len(nrow(x)))
    colnames(x) <- paste0('cell',seq_len(ncol(x)))
    o <- CreateSeuratObject(x)
    o$patient <- rep(c('synthetic_a','synthetic_b'),length.out=ncol(o))
    o$lineage <- factor(rep(c('CD8','CD4','NK'),length.out=ncol(o)))
    o$neutral_flag <- rep(c(TRUE,FALSE,NA),length.out=ncol(o))
    qs2::qs_save(o,file.path(d,'native.qs2'))
    qs2::qs_to_rds(file.path(d,'native.qs2'),file.path(d,'native.rds'),compress_level=0)
    saveRDS(o,file.path(d,'xdr.rds'))
    write.csv(o@meta.data,file.path(d,'expected.csv'),row.names=TRUE,na='')
    '''
    script = directory/'fixture.R'; script.write_text(r_source)
    subprocess.run(['Rscript',str(script),str(directory)],cwd=ROOT,check=True)
    expected = pd.read_csv(directory/'expected.csv',index_col=0)
    expected['neutral_flag'] = expected['neutral_flag'].astype('boolean')
    checks = {}
    for encoding in ['native','xdr']:
        actual,audit = extract(directory/f'{encoding}.rds',max_elements=1000)
        pd.testing.assert_frame_equal(actual,expected,check_dtype=False,check_categorical=False)
        assert len(audit['discarded_numeric_arrays']) >= 2
        try:
            extract(directory/f'{encoding}.rds',max_elements=250)
        except AssertionError as error:
            assert 'Metadata contains a discarded array' in str(error)
        else:
            raise AssertionError('Unsafe metadata truncation was not rejected')
        checks[encoding] = {'all_metadata_values_match_R':True,'unsafe_limit_rejected':True,
                            'discarded_numeric_arrays':len(audit['discarded_numeric_arrays'])}
    (ROOT/'reports/evidence/cohort-expansion/metadata-stream-validation.json').write_text(
        json.dumps({'synthetic_only':True,'rows':len(expected),'checks':checks},indent=2)+'\n')
    print('PASS: both R serialization formats preserve every metadata value; unsafe truncation rejected.')
