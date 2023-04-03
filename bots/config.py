import tweepy as tw
import datetime



# tweepy-bots/bots/config.py
import tweepy
import logging
import os

logger = logging.getLogger()

def create_api():
    consumer_key = os.getenv("uxzUXjkFEbR883ShiKPtKt58x")
    consumer_secret = os.getenv("6F38WQtUfvflFSPnt7mvqFybgiBK53T2WuYrz0Ug6VdSgbDJYP")
    access_token = os.getenv("1340153593526185985-lnVeI8a9L8c301NKBVWSzMhkvbp7GL")
    access_token_secret = os.getenv("QPFL731pUT7cQgwaeqsScuMzbghHzxCPaFQtJUKO3orj3")

    auth = tweepy.OAuthHandler(consumer_key, consumer_secret)
    auth.set_access_token(access_token, access_token_secret)
    api = tweepy.API(auth, wait_on_rate_limit=True, 
        wait_on_rate_limit_notify=True)
    try:
        api.verify_credentials()
    except Exception as e:
        logger.error("Error creating API", exc_info=True)
        raise e
    logger.info("API created")
    return api

