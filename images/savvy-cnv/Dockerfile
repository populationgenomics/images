FROM broadinstitute/gatk:latest

ARG VERSION=${VERSION}

#COPY README* *build.d /build.d/
#RUN bash -c 'for i in $(ls /build.d/*.sh 2>/dev/null | sort) ; do source $i ; done'

RUN apt-get update

RUN apt-get install -y git

WORKDIR /home/gradle/

RUN git clone https://github.com/rdemolgen/SavvySuite.git \
&& cd SavvySuite && javac *.java

ENV CLASSPATH="/home/gradle/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar:/home/gradle/SavvySuite:${CLASSPATH}"

CMD []
