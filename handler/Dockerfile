FROM python:3.11-slim

RUN mkdir /app
WORKDIR /app

# Install necessary Python packages
RUN python -m pip install --upgrade pip
RUN pip install requests boto3

COPY PAHandler.py /app/PAHandler.py

# Command to run the Python script
CMD ["python", "/app/QCHandler.py"]