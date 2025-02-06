venv:
	python3.12 -m venv .venv

r:
	./.venv/bin/pip install -r requirements.txt

run:
	./.venv/bin/python3.12 main.py