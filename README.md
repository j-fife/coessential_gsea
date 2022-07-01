# Coessential GSEA    

 
 
### clone the repositort
             
```bash
git clone https://github.com/j-fife/coessential_gsea.git .
cd ./src/
```

### Usage

```python 
from gsea import CoessentialGSEA

input_file = "example_data.csv"
c = CoessentialGSEA(input_file)
c.run(d = 0.2, permutation_num = 10000)
```

