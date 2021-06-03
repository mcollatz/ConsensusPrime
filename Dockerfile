FROM docker.io/matthewfeickert/docker-python3-ubuntu:latest


RUN /bin/bash


RUN sudo apt-get update -y
RUN sudo apt-get install -y mafft
RUN sudo apt-get install -y primer3
RUN pip3.8 install pandas

ADD consensus_prime.py /consensus_prime.py
CMD [ "/consensus_prime.py" ,"-h" ]
ENTRYPOINT [ "python3.8" ]

