clustalo-python
===============

This is just a simple Python wrapper around Clustal Omega
(http://www.clustal.org/omega/), used internally at Benchling but casually open
source, if it helps anybody. Also available via::

  pip install clustalo

Requires `libclustalo` installed, version 1.2.0. On Linux, it's recommended you
build and install it yourself::

  cd clustal-omega-1.2.0
  ./configure --with-pic --with-openmp
  make && sudo make install

before installing this package.

Support for OSX is not quite there yet, mainly because OpenMP isn't supported
on clang. Homebrew is still on 1.1.0, so you'll need to compile and install
clustalo 1.2.0 yourself (`--without-openmp`). You'll also need to set
`OPENMP_DISABLED=true` in env vars before running build/install.

Usage
-----
::

  from clustalo import clustalo
  input = {
      'seq1': 'AAATCGGAAA',
      'seq2': 'CGGA'
  }
  aligned = clustalo(input)
  # aligned is a dict of aligned sequences:
  #   seq1: AAATCGGAAA
  #   seq2: ----CGGA--

At the moment, input sequences are assumed to not be aligned (i.e. there is no
dealign option). See ``clustalo.clustalo.__doc__`` or file ``clustaslo/clustalo.c``
for documentation.

