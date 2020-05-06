FROM python:3.8.2

RUN apt-get update -qqy \
    && apt-get install -qqy --no-install-recommends \
        libargtable2-dev \
        python-dev

RUN wget -q http://www.clustal.org/omega/clustal-omega-1.2.4.tar.gz \
    && tar -xzvf clustal-omega-1.2.4.tar.gz \
    && cd clustal-omega-1.2.4 \
    && ./configure --with-pic --with-openmp \
    && make \
    && make install

RUN python3 -m pip install --user --upgrade \
        pip \
        setuptools \
        wheel

RUN wget -q https://bootstrap.pypa.io/get-pip.py \
    && python2 get-pip.py --user \
    && python2 -m pip install --user --upgrade \
        pip \
        setuptools \
        wheel

ADD . /src

RUN cd /src \
    && python3 setup.py sdist bdist_wheel \
    && python2 setup.py sdist bdist_wheel
