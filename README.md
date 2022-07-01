# Coessential GSEA    

 
 
### Clone the repositort
             
```bash
git clone https://github.com/j-fife/coessential_gsea.git .
pip install -r requirements.txt
cd ./src/
```



### Usage

```python 
from gsea_tools import CoessentialGSEA

input_file = "example_data.csv"
c = CoessentialGSEA(input_file)
c.run(d = 0.2, permutation_num = 10000)
```

