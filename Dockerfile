FROM ensemblorg/ensembl-vep:release_113.3

USER root
RUN apt-get update  && apt-get install -y less
RUN apt install -y python3-pip wget git

RUN mkdir /opt/vep/ronavep

COPY requirements.txt /opt/vep/ronavep

WORKDIR /opt/vep/ronavep

RUN pip3 install setuptools==58  #Downgrading setuptools to install pyvcf
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

RUN pip3 install git+https://github.com/andersen-lab/Freyja@2ea438c590b33e4c105e422f0c03f4a2e84d7870

CMD python3 /opt/vep/ronavep/main.py