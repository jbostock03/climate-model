FROM jupyter/base-notebook:python-3.10

# Install Python packages
COPY .binder/requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Copy repo contents
COPY . /home/jovyan/

# Make all notebooks trusted
RUN find /home/jovyan -name "*.ipynb" -exec jupyter trust {} \;

# Make custom start script executable
COPY .binder/start /usr/local/bin/start
RUN chmod +x /usr/local/bin/start

# Use it as the entrypoint if needed
ENTRYPOINT ["/usr/local/bin/start"]

# Default command (launch JupyterLab)
CMD ["start.sh"]
