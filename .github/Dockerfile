FROM python:3.6

COPY setup.py ./

RUN pip install --upgrade pip && \
    pip install --no-cache-dir

COPY . .

CMD ["bash", "genomepy"]
