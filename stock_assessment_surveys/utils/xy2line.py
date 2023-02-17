import os
import shutil

import fiona

def xy2line(df):
    # get a distinct list of years
    years = sorted(df["year__annee"].unique().tolist())

    schema = {
        'geometry': 'LineString',
        'properties': [('YearAnnee', 'str'),('StartDepart', 'str'),('EndFin', 'str')]
    }

    wms_dir = "./WMS"
    shapefile = "CoverageByYear.shp"
    shapefile_dir = os.path.join(wms_dir, shapefile)
    with fiona.open(shapefile_dir, mode='w', driver='ESRI Shapefile', schema=schema, crs="EPSG:4326") as lineShp:
        for year in years:
            mask = df["year__annee"] == year
            yearly_df = df[mask]
            coords = list()
            for index, row in yearly_df.iterrows():
                coords.append((row["longitude"], row["latitude"]))
            start = yearly_df["datetime"].min()
            end = yearly_df["datetime"].max()

            rowDict = {
                'geometry': {'type': 'LineString',
                             'coordinates': coords},
                'properties': {'YearAnnee': year, 'StartDepart': start.strftime("%Y-%m-%d"), 'EndFin': end.strftime("%Y-%m-%d")},
            }
            lineShp.write(rowDict)

    # add all files into a zip
    shutil.make_archive(shapefile, 'zip', "./WMS")