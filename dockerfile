# FROM python:3.7.6
FROM python:buster

# RUN mkdir /home ; exit 0
RUN mkdir /app
WORKDIR /app
COPY requirements.txt .
RUN python -m pip install --upgrade pip
RUN pip install -r requirements.txt
RUN ls 
COPY . .
EXPOSE 8080
CMD ["python", "app.py"]