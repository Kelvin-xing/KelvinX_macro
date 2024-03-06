Here's a Python code that reads the given series of strings, classifies entity name, alias, subordinate/affiliate, and address, and concludes patterns of each classification in regular expressions:

```python
import re

# Regular expressions for entity name, alias, subordinate/affiliate, and address
name_pattern = r'([A-Z][a-z]+(?: [A-Z][a-z]+)*(?: \([A-Z]+\))?(?: [A-Z][a-z]+)*)'
alias_pattern = r'(?:a\.k\.a\.|aka|also known as|alias|the following two aliases:)\s*(.*?)(?=(?:,|;|\n))'
subordinate_pattern = r'(?:,|and)\s*(.*?)(?=(?:,|;|\n))'
address_pattern = r'([A-Z][a-z]+(?: [A-Z][a-z]+)*(?:, [A-Z][a-z]+)*)'

# Loop over each element of the series
for string in series:
    # Extract entity name
    entity_name = re.search(name_pattern, string).group(1)
    
    # Extract aliases
    aliases = re.findall(alias_pattern, string)
    
    # Extract subordinates/affiliates
    subordinates = re.findall(subordinate_pattern, string)
    
    # Extract address
    addresses = re.findall(address_pattern, string)
    
    # Print the results
    print("Entity Name:", entity_name)
    print("Aliases:", aliases)
    print("Subordinates/Affiliates:", subordinates)
    print("Addresses:", addresses)
    print("\n")
```

This code uses regular expressions to extract the entity name, aliases, subordinates/affiliates, and addresses from each string in the series. It then prints the results for each string.