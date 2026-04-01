.PHONY: install test run benchmark plots clean info

install:
	bash setup.sh

test:
	python -m pytest tests/ -v

run: test

benchmark:
	python -m benchmarks.bench_runner

plots:
	python -m benchmarks.plot_results

info:
	python -m benchmarks.hardware_info

clean:
	rm -rf results/*.json results/*.png
	rm -rf __pycache__ genec/__pycache__ tests/__pycache__ benchmarks/__pycache__
	rm -rf .pytest_cache tests/.pytest_cache
