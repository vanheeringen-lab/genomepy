FROM python:3.6

COPY setup.py ./

RUN python setup.py install

COPY . .

CMD ["bash", "genomepy"]
