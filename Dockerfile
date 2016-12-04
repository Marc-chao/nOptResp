FROM python:2.7

WORKDIR /nOptResp

RUN pip install -r requirements.txt

CMD python -m unittest discover tests
CMD python exec/solve_E0dep.py
