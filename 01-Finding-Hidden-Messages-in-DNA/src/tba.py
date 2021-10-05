
class tba:

   def __init__(self) -> None:
       print("tba constructor")

   def pattern_count(self, text: str, pattern: str) -> int:

      count = 0
      pattern_len = len(pattern)
      for i in range(0, len(text) - len(pattern)):
         j = i + pattern_len
         if text[i:j] == pattern:
            count +=1

      return count



