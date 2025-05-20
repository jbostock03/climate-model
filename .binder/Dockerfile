FROM jupyter/base-notebook:python-3.10

# Install Python packages
COPY .binder/requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Copy repo contents
COPY . /home/jovyan/

# Make all notebooks trusted
RUN find /home/jovyan -name "*.ipynb" -exec jupyter trust {} \;

# Execute all notebooks
RUN find /home/jovyan -name "*.ipynb" -exec jupyter nbconvert --to notebook --execute --inplace --allow-errors {} \;

# Default command (launch JupyterLab)
CMD ["start.sh"]
