# meteokav-ecm9

ECMWF IFS 9 km (0.1°) total precipitation over Greece, as 3-hour totals for 96 h, for the precipitation map of
[meteokav.gr](https://www.meteokav.gr).

- `ecm9/build.py`: downloads only the precipitation field of each forecast step from ECMWF's open data (using the
  `.index` byte ranges), decodes it with eccodes, crops the Greece box and writes `latest.json` plus `p<run>.png`.
  The comments in the file describe the format.
- `.github/workflows/ecm9.yml`: runs every 30 minutes. When a new run (00/06/12/18 UTC) is out, it builds and
  force-pushes the result to the `data` branch.
- The website pulls `data/latest.json` and the picture from `raw.githubusercontent.com`.

Data: © ECMWF, open data under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Source: https://www.ecmwf.int/en/forecasts/datasets/open-data
