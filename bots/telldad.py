import tweepy as tw
import datetime

# Authenticate to Twitter
auth = tw.OAuthHandler("uxzUXjkFEbR883ShiKPtKt58x",
	"6F38WQtUfvflFSPnt7mvqFybgiBK53T2WuYrz0Ug6VdSgbDJYP")

auth.set_access_token("1340153593526185985-lnVeI8a9L8c301NKBVWSzMhkvbp7GL",
	"QPFL731pUT7cQgwaeqsScuMzbghHzxCPaFQtJUKO3orj3")

api = tw.API(auth)

try:
   api.verify_credentials()
   print("Authentication OK")
except:
   print("Error during authentication")


#public tweet test==================================================================
def publictweet_1():
    if datetime.date.today().weekday() == 0:
        tweettopublish = 'Hi everyone, today is Monday.   #Monday '
    if datetime.date.today().weekday() == 1:
        tweettopublish = 'Enjoy your Tuesday.  #Tuesday'
    if datetime.date.today().weekday() == 2:
        tweettopublish = 'Third week of the Week. #Wednesday'
    if datetime.date.today().weekday() == 3:
        tweettopublish = 'Thursday. I cannot wait for the Weekend'
    if datetime.date.today().weekday() == 4:
        tweettopublish = 'Friday...Finally'
    if datetime.date.today().weekday() == 5:
        tweettopublish = 'Great it is Saturday #weekend #Saturday'
    if datetime.date.today().weekday() == 6:
        tweettopublish = 'Sunday morning...#Weekend #enjoy '

    api.update_status(tweettopublish)
    print(tweettopublish)

print(datetime.date.today().weekday())

#====================================================================================
def publictweet_2():
  screen_name = "Can't wait to tell my dad"
  c = tw.Cursor(api.followers, screen_name)
  count = 0
  for follower in c.items():
   count += 1
  print("I can't wait to tell my dad that I have " + str(count) + " followers on Twitter")
  twe = "I can't wait to tell my dad that I have " + str(count) + " followers on Twitter"

  api.update_status(twe)






# execute

#publictweet_1()
publictweet_2()
