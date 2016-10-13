import multiprocessing

def a(l):
    return l+1
tasks = []
for i in range(4):
    tasks.append((i, ))

pool = multiprocessing.Pool(2)
results = [pool.apply_async(a, t) for t in tasks]
results
for result in results:
    b = result.get()
    print b