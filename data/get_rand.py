import random

#for i in range(2*520*520*520):
#    print(random.random()*100)
random.seed(233)
cnt = 10
for i in range(cnt):
    print("%.2f"%(random.random()*100))
