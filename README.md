clustalo-python
===============

This is just a simiple Python wrapper around Clustal Omega
(http://www.clustal.org/omega/), used internally at Benchling but casually open
source if it helps anybody. Also available via

```
pip install clustalo
```

To use,

```
from clustalo import clustalo
input = {
    'seq1': 'AAATCGGAAA',
    'seq2': 'CGGA'
}
aligned = clustalo(input)
# aligned is a dict of aligned sequences:
#   seq1: AAATCGGAAA
#   seq2: ----CGGA--
```

At the moment, input sequences are assumed to not be aligned (i.e. there is no
dealign option). The following optional keyword parameters are supported:

- **seqtype**: clustalo.DNA, clustalo.RNA, clustalo.PROTEIN, or clustalo.UNKNOWN (to have it guess for you)
- **mbed_guide_tree**: True, False (default True)
- **mbed_iteration**: True, False (default True)
- **num_combined_iterations**: int (default 0)
- **max_guidetree_iterations**: int (default INT_MAX)
- **max_hmm_iterations**: int (default INT_MAX)


