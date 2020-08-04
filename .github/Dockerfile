FROM python:3.6

COPY setup.py ./
COPY genomepy/__about__.py ./
COPY README.md ./

RUN python setup.py install

COPY . .

CMD ["bash", "genomepy"]
