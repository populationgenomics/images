FROM python:3.8 
RUN pip install cmake && pip install make

FROM ubuntu
RUN mkdir build; cd build
RUN cmake ..; make