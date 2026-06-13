.PHONY: setup pipeline dashboard

setup:
	python -m pip install -r requirements.txt

pipeline:
	python load_data.py
	python analysis.py

dashboard:
	python -m streamlit run dashboard.py
