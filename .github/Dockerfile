FROM python:3.6

COPY setup.py ./

RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r setup.py

COPY . .

CMD ["bash", "genomepy"]
