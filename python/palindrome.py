#  A palindromic number reads the same both ways. The largest palindrome made
#  from the product of two 2-digit numbers is 9009 = 91 99.

#  Find the largest palindrome made from the product of two 3-digit numbers.

def find_palindrome(number):
    numList = map(int, str(number))
    length = len(numList)
    i = 1
    j = 0
    k = -1
    palindrome = True
    while i <= len(numList)/2 and palindrome == True:
        if numList[j] != numList[k]:
            palindrome = False
        i += 1
        j += 1
        k -= 1
    return [number, palindrome]

i = 999
j = 999
palindrome = False
nums = []
while i > 99 and j > 99: # and palindrome == False:
    number = i * j
    palCheck = find_palindrome(number)
    if palCheck[1] == True:
        nums.append(palCheck[0])
    if i < 101:
        i = 999
        j -= 1
    else:
        i -= 1

print max(nums)
        
