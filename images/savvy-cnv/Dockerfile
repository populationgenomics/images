from gatk

RUN git clone https://github.com/rdemolgen/SavvySuite.git

ENV CLASSPATH="/home/gradle/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar:/home/gradle/SavvySuite:${CLASSPATH}"

WORKDIR /home/gradle/SavvySuite

RUN javac *.java

CMD []
