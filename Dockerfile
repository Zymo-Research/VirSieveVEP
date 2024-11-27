FROM ensemblorg/ensembl-vep:release_103

USER root
RUN apt-get update  && apt-get install -y less
RUN apt install -y python3-pip
RUN ln -s /usr/bin/python3 /usr/bin/python

RUN mkdir /opt/vep/ronavep

COPY requirements.txt /opt/vep/ronavep

WORKDIR /opt/vep/ronavep

RUN pip3 install -r requirements.txt

COPY references /opt/vep/ronavep/references

WORKDIR /opt/vep/ronavep/references

RUN /bin/bash prepRef.sh

COPY ./*.py /opt/vep/ronavep/

COPY cvaSupport /opt/vep/ronavep/cvaSupport

RUN chown -R vep /opt/vep/ronavep

USER vep

ENV PYTHONUNBUFFERED=1

WORKDIR /opt/vep

CMD python3 /opt/vep/ronavep/main.py