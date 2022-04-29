
## Example 1
##  an error will occur because
##    python has a maximum recursion depth of 1000

# def hello():
#   print('Hello, world!')
#    hello()
# hello()


## Example 2 

def hello(count):
    if count == 0: # ending condition
        return

    print('Hello, world!', count)

    count -= 1 

    hello(count) 

hello(5)
