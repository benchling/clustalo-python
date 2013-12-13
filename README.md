# clustalo-python

This is just a simiple Python wrapper around Clustal Omega
(http://www.clustal.org/omega/), used internally at Benchling but casually open
source if it helps anybody. Also available via

```
pip install clustalo
```

Requires `libclustalo` installed, version 1.2.0. On Linux, it's recommended you build and
install it yourself:

```
cd clustal-omega-1.2.0
./configure --with-pic --with-openmp
make && sudo make install
```

Homebrew (might) be still on 1.1.0, so if on Mac you might need to build it
yourself also. OpenMP isn't available for clang, so you can either install gcc:

```
brew tap homebrew/versions
brew install gcc48
```

or just install without it. Just know that `num_threads` won't do anything.

## Usage

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

- **seqtype**: `clustalo.DNA`, `clustalo.RNA`, `clustalo.PROTEIN`, or `clustalo.UNKNOWN` (to have it guess for you)
- **mbed_guide_tree**: `True`, `False` (default `True`)
- **mbed_iteration**: `True`, `False` (default `True`)
- **num_combined_iterations**: `int` (default `0`)
- **max_guidetree_iterations**: `int` (default `INT_MAX`)
- **max_hmm_iterations**: `int` (default `INT_MAX`)
- **num_threads**: `int` (default `1`) - needs OpenMP to do anything

