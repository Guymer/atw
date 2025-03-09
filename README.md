# Around The World (ATW)

!["mypy" GitHub Action Status](https://github.com/Guymer/atw/actions/workflows/mypy.yaml/badge.svg) !["pylint" GitHub Action Status](https://github.com/Guymer/atw/actions/workflows/pylint.yaml/badge.svg)

This project aims to show an interesting route for you to circumnavigate planet Earth.

## Workflow

1. Find the most interesting great circle to your antipode (by running [step1_pickGreatCircle.py](step1_pickGreatCircle.py))
2. Show all the countries that your chosen great circle intersects (by running [step2_showGreatCircle.py](step2_showGreatCircle.py))
3. Animate your journey (by running [step3_animateRoute.py](step3_animateRoute.py))

## Dependencies

ATW requires the following Python modules to be installed and available in your `PYTHONPATH`.

* [cartopy](https://pypi.org/project/Cartopy/)
* [geojson](https://pypi.org/project/geojson/)
* [matplotlib](https://pypi.org/project/matplotlib/)
* [numpy](https://pypi.org/project/numpy/)
* [pyguymer3](https://github.com/Guymer/PyGuymer3)
* [shapely](https://pypi.org/project/Shapely/)
