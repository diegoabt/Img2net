[![CodeFactor](https://www.codefactor.io/repository/github/diegoabt/img2net/badge)](https://www.codefactor.io/repository/github/diegoabt/img2net)

# Img2net: from images to networks

Python implementation of the _img2net_ algorithm proposed in:

- [1]      Baptista Diego and De Bacco Caterina. 2021. <span style="color:blue">[Principled network extraction from images](https://doi.org/10.1098/rsos.210025) *R. Soc. open sci.* 8210025210025 .

<img src="figures/pre_extracted__angolan_river_crop.png" width="40%" alt="hi" class="inline"/>

If you use this tool please cite [1].

## What's included?

- `code/` :
    + `img2net-core/`: main files used to extract the networks from the images.
    + `metrics/`: scripts used to compute the metrics exhibited in [1].
    + `test/`: test files.
- `figures/` : README figures.


## Requirements

This repository uses some of the functions proposed on the repository _Nextrout_[2], thus, it is necessary to clone it:

- [Nextrout](https://github.com/Danielaleite/Nextrout)

That repository itself needs some dependecies coming from some GitLab repositories, so you will need an [account](https://about.gitlab.com/).

Once _Nextrout_ is cloned, the next step is to install the other python packages found in `requirements.txt`:

```
pip install -r requirements.txt
```

Please locate  `Img2net` and `Nextrout` at the same level, i.e.,

- folder
    + ...
    + Image2Net
    + Nextrout
    + ...

## Test

To test the tool, please go the folder `test` inside `code` folder and execute

```bash
python test_img2net.py
```
This will execute _img2net_ on a test image contained in the `input` folder; results will be stored in the folder `runs`.

## References

- [2] Baptista, D., Leite, D., Facca, E. et al. *Network extraction by routing optimization. Sci Rep 10, 20806 (2020)*. [https://doi.org/10.1038/s41598-020-77064-4](https://doi.org/10.1038/s41598-020-77064-4)
