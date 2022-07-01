### Coessential GSEA        
                                                               
             
```bash
git clone https://github.com/j-fife/coessential_gsea.git .
cd ./src/
```

```python 
from gsea import CoessentialGSEA

input_file = "example_data.csv"
c = CoessentialGSEA(input_file)
c.run(d = 0.2, permutation_num = 10000)
```
