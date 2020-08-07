# BASE IMAGE
FROM biocontainers/biocontainers:latest
# METADATA
# LABELs ...
# Reference to developers (Kim)
# MAINTAINER
MAINTAINER VIB Bioinformatics Core (tuur.muyldermans@vib.be)
USER biodocker

ADD . /data

ENV DST=/home/biodocker/bin
ENV WRKDR=/data

# WORKDIR /data

# RUN wget $URL/$ZIP -O $DST/$NAME \
#     && tar -xzvf $DST/$NAME -C $DST/
    # mv $DST/$FOLDER $WRKDR  
RUN conda create -y -c bioconda -n tmp-env python=2.7 numpy=1.16.5 bamtools=2.5.1
#SHELL ["source", "activate", "tmp-env"]
#SHELL ["conda", "run", "-n", "tmp-env", "/bin/bash", "-c"] 
#RUN conda create -y -c bioconda -n tmp-env python=2.7 numpy=1.16.5 bamtools=2.5.1 && source activate tmp-env 
RUN echo "source activate tmp-env" > ~/.bashrc
ENV PATH /opt/conda/envs/tmp-env:$PATH
RUN cp -rf $DST/$FOLDER/* $WRKDR
#ENTRYPOINT ["source", "activate", "tmp-env", "python", "/home/biodocker/bin/RDDpred-1.1.0/RDDpred.py" ]
#ENTRYPOINT ["python", "/home/biodocker/bin/RDDpred-1.1.0/RDDpred.py"]
#ENV PATH /home/biodocker/bin/$FOLDER:$PATH