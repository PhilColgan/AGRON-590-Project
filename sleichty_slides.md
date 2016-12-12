sleichty_slides
========================================================
author: 
date: 
autosize: true

Data Dictionary Format
========================================================
Sections for: 
- column name in spreadsheet
- full column name 
  - ex. veg_class = Vegetation Class
- descriptions of process
  - ex. Gas measurements in micromoles per square meter
Packages Used for Data Dictionary
========================================================

- Kable
  - Creates rectangular table
  - kable(x, format, digits = getOption("digits"), row.names = NA, col.names = NA, align, 
    caption = NULL, format.args = list(), escape = TRUE, ...)
  - knitr::kable(head(data_dictionary))
Special Changes for Output
=======================================================

```
[1] "$\\beta$-N-acetylglucosaminidase Activity"
```
"$\\beta$-N-acetylglucosaminidase Activity"
  - Code to produce Beta symbol
  - Important for correct enzyme assignment
Data Dictionary Output
========================================================

|name         |plot_name    |description                    |
|:------------|:------------|:------------------------------|
|plot         |Plot         |Plot sample was taken from     |
|treatment    |Treatment    |Treatment number plot received |
|sample_month |Sample Month |Month the sample was taken     |
|sample_year  |Sample Year  |Year the sample was taken      |
|crop         |Crop         |Cropping system                |
|sample_block |Sample Block |Experimental block             |


