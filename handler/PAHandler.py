import requests
import boto3
import json
import os

client = boto3.client('sqs', region_name='us-west-2')
sns = boto3.client('sns', region_name='us-west-2')
queue_url = os.environ.get("SQS_QUEUE_URL")
env = os.environ.get("ENVIRONMENT", "dev")
print(f"Cellborg Troop - PA Python container running in {env} environment")

if env == "dev":
    SNS_TOPIC = 'arn:aws:sns:us-west-2:865984939637:QCCompleteTopic' #there is no dev environment for PA
else:
    SNS_TOPIC = f'arn:aws:sns:us-west-2:865984939637:PAComplete-{env}-Topic'

def send_sns(data):
    topic_arn = SNS_TOPIC
    print(f'Sending {data} SNS message to {SNS_TOPIC}')
    response = sns.publish(
        TopicArn = topic_arn,
        Message = json.dumps(data),
    )
    return response

def send_request(endpoint, data):
    url = f"http://127.0.0.1:8001{endpoint}"
    response = requests.post(url, json = data, headers = {'Content-Type': 'application/json'})
    return response.json()

def send_shutdown_request(user=None, project=None):
    url = f"http://127.0.0.1:8001/shutdown"
    try:
        if user and project:
            response = requests.post(url, json = {"user": user, "project": project}, headers = {'Content-Type': 'application/json'})
        else:
            response = requests.post(url, json = {"status": "complete"}, headers = {'Content-Type': 'application/json'})
        response.raise_for_status()
    except requests.ConnectionError:
        print("Connection was closed by the server (expected behavior during shutdown).")
    except requests.RequestException as e:
        print(f"An error occurred: {e}")
    
MAX_COUNT = 1000
currentCount=0

#send notification that task is running to api
while True:
    response = client.receive_message(QueueUrl=queue_url, MaxNumberOfMessages=10, WaitTimeSeconds=10, VisibilityTimeout=900)
    print("queueurl=",queue_url)
    print(response)

    if currentCount>= MAX_COUNT:
        print("Server hashit timeout, shutting down...")
        send_shutdown_request()

    if 'Messages' not in response:
        print("No Message in ",queue_url, "topic:",SNS_TOPIC)
        currentCount+=1
        if currentCount%10==0:
            print(currentCount)
        continue

    for message in response['Messages']:
        try:
            queen_service_request_raw_data = message['Body']
            queen_service_request = json.loads(queen_service_request_raw_data)
            print(queen_service_request)
            request_type = queen_service_request["requestType"]
            project = queen_service_request["project"]
            user = queen_service_request["user"]


            if request_type == "initializeProject":
                print("Initializing annData object")
                datasets = queen_service_request["datasets"]

                response = send_request('/init_endpoint', {"user": user, "project": project, "datasets": datasets})
                if response["success"]:
                    print("Initializing Project Successful... Sending SNS message")
                    project_initialized = True
                    data = {
                        "user": user, 
                        "project": project, 
                        "stage": "initialized"
                    }
                    response = send_sns(data)
                    print(response)
            elif request_type == "clustering":
                print("Beginning clustering now...")
                resolution = queen_service_request['resolution']

                response = send_request('/clustering', {"user":user, "project":project, "resolution":resolution})
                if response['success']:
                    print("Clustering was successful...")
                    data = {
                        "user":user,
                        "project":project,
                        "clusters": response["clusters"],
                        "stage":"cluster"
                    }
                    response = send_sns(data)
                    print(response)
                    
            elif request_type == "gene_expression":
                print("Beginning gene expression now...")
                gene_list = queen_service_request['gene_list']
                response = send_request("/gene_expression", {"user":user, "project":project, "gene_list":gene_list})
                if response['success']:
                    print("gene expression completed")
                    data = {
                        "user":user,
                        "project":project,
                        "stage": "gene_expression"
                    }
                    response= send_sns(data)
                    print(response)

            elif request_type == "annotations":
                print("Beginning annotations now...")
                anno = queen_service_request['annotations']
                response = send_request("/annotations", {"user":user, "project":project, "annotations":anno})
                if response['success']:
                    print("Annotations happened successfully")
                    data = {
                        "user":user,
                        "project":project,
                        "stage":"annotations"
                    }
                
                    response = send_sns(data)
                    print(response)
            elif request_type == "killServer":

                print("Saving adata and shutting down the PA server...")
                send_shutdown_request(user, project)
                print("Shutting down the python handler")
                exit(0)
            
        except Exception as e:
            print("Error:", str(e))
        finally:
            client.delete_message(QueueUrl=queue_url, ReceiptHandle=message['ReceiptHandle'])