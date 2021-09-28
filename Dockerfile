FROM ubuntu:18.04

RUN mkdir -p /opt/sent2vec_api/  # make WORKDIR

RUN apt update
RUN apt install -y libevent-pthreads-2.1-6
RUN apt-get update -y
RUN apt-get install -y python3.8  # install python 3.8
RUN apt-get install -y python3-pip python3.8-dev build-essential  # install python build tool
RUN apt-get update -y
RUN apt-get install -y git
RUN ln -s /usr/bin/python3.8 /usr/bin/python  # Setup python 3.8 as default

# Install Sent2Vec
WORKDIR /opt/sent2vec_api
RUN python3.8 -m pip install --upgrade pip
RUN git clone https://github.com/epfml/sent2vec.git
WORKDIR sent2vec
RUN make
RUN python3.8 -m pip install Cython
RUN python3.8 -m pip install numpy
RUN python3.8 -m pip install .
WORKDIR /opt/sent2vec_api

# Prepare Sent2Vec pipline
ADD sent2vec_pipeline.py /opt/sent2vec_api

# Add all reference scripts
ADD ../fetch_pmid_title_abstract.py /opt/sent2vec_api
ADD ../hmac_signature.py /opt/sent2vec_api