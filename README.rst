clustalo-python
===============

This is just a simple Python wrapper around Clustal Omega
(http://www.clustal.org/omega/), used internally at Benchling but casually open
source, if it helps anybody. Also available via::

  pip install clustalo

Requires `libclustalo` installed, version 1.2.0. On Linux, it's recommended you build and
install it yourself::

  cd clustal-omega-1.2.0
  ./configure --with-pic --with-openmp
  make && sudo make install

Homebrew (might) be still on 1.1.0, so if on Mac you might need to build it
yourself also. OpenMP isn't available for clang, so you can either install gcc::

  brew tap homebrew/versions
  brew install gcc48

or just install without it. Just know that `num_threads` won't do anything.

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

