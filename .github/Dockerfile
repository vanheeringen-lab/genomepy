FROM python:3.6

COPY . ./

RUN python setup.py install

ENTRYPOINT ["genomepy"]
CMD ["--help"]
