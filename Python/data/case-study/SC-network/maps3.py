import os   # Import the os module: provides a way of using operating system dependent functionality
import osmnx as ox  # Import the osmnx module: a Python package that lets you download spatial data from OpenStreetMap
import geopandas as gpd # Import the geopandas module: a Python package that lets you work with geospatial data
import matplotlib.pyplot as plt # Import the matplotlib module: a Python 2D plotting library
import folium
import math


def get_map(place_name):
    # let's first check if the data is already downloaded
    if not os.path.exists(f"{place_name}.graphml"): # If the data is not downloaded
        G = ox.graph_from_place(place_name, network_type='all') # Download the data
        ox.save_graphml(G, f"{place_name}.graphml") # Save the data to a file
    else:
        G = ox.load_graphml(f"{place_name}.graphml")


    #let's check if the boundaries are already downloaded
    if not os.path.exists(f"{place_name}.geojson"):
        gdf = ox.geocode_to_gdf(place_name) # Download the boundaries of the specified place
        gdf.to_file(f"{place_name}.geojson", driver='GeoJSON') # Save the boundaries to a file
    else:
        gdf = gpd.read_file(f"{place_name}.geojson")

    return G, gdf
    

# Define the place you want to download the map for
place_name = "South Carolina, USA"
G, gdf = get_map(place_name)
# coastal counties in South Carolina
counties = ['Horry County', 'Georgetown County', 'Charleston County', 'Beaufort County', 'Jasper County', 'Colleton County', 'Dorchester County', 'Berkeley County']








def get_data(save_data=False):
    # let's create a DataFrame to store the airports
    airports_df          = gpd.GeoDataFrame()
    schools_df           = gpd.GeoDataFrame()
    religious_centers_df = gpd.GeoDataFrame()
    for county in counties:
        tags = {'aeroway': True}
        county_name = f"{county}, South Carolina, USA"
        airports = ox.features_from_place(county_name, tags)
        # let's filter out only the points (i.e., row.geometry.geom_type == 'Point') that are airports (aerodromes), and drop any rows with missing names
        
        airports = airports[airports.apply(lambda row: 'geometry' in row and row.geometry.geom_type == 'Point', axis=1)]
        airports = airports[airports['aeroway'] == 'aerodrome'] 
        airports = airports.dropna(subset=['name'])
        airports['county'] = county
        airports_df = airports_df.append(airports)
        
        tags = {'amenity': 'school'}
        county_name = f"{county}, South Carolina, USA"
        schools = ox.features_from_place(county_name, tags)
        # let's filter out only the points (i.e., row.geometry.geom_type == 'Point') that are schools, and drop any rows with missing names
        schools = schools[schools.apply(lambda row: 'geometry' in row and row.geometry.geom_type == 'Point', axis=1)]
        schools = schools.dropna(subset=['name'])
        schools['county'] = county
        schools_df = schools_df.append(schools)
        
        
        tags = {'amenity': 'place_of_worship'}
        county_name = f"{county}, South Carolina, USA"
        religious_centers = ox.features_from_place(county_name, tags)
        # let's filter out only the points (i.e., row.geometry.geom_type == 'Point') that are religious centers, and drop any rows with missing names
        religious_centers = religious_centers[religious_centers.apply(lambda row: 'geometry' in row and row.geometry.geom_type == 'Point', axis=1)]
        religious_centers = religious_centers.dropna(subset=['name'])
        religious_centers['county'] = county
        religious_centers_df = religious_centers_df.append(religious_centers)
                
        
    # let's save the data to csv files
    if save_data:
        airports_df.to_csv('airports.csv')
        schools_df.to_csv('schools.csv');
        religious_centers_df.to_csv('religious_centers.csv')
    
    
    return airports_df, schools_df, religious_centers_df

#airports_df, schools_df, religious_centers_df = get_data(save_data=True);


# function to add a county to the map with a blue highlight and a marker
def add_county_to_map(county, map_obj, show_airports=True, show_schools=True, show_religious_centers=False):
    county_name = f"{county}, South Carolina, USA"
    cgdf = ox.geocode_to_gdf(county_name)
    
    # Add the county boundary with blue highlight
    folium.GeoJson(
        cgdf,
        style_function=lambda feature: {
            'fillColor': 'darkblue',
            'color': 'darkblue',
            'weight': 2,
            'fillOpacity': 0.3
        }
    ).add_to(map_obj)
    

    if not (show_airports or show_schools or show_religious_centers):
        cgdf.geometry = cgdf.geometry.to_crs('EPSG:3857')
        centroid = cgdf.geometry.centroid
        centroid = centroid.to_crs('EPSG:4326')
        folium.Marker(
            location=[centroid.y, centroid.x],
            popup=county,
            tooltip=county,
            icon=folium.Icon(color='darkblue', icon='info-sign')
        ).add_to(map_obj)

    if show_airports:
        # Let's mark the airports in the county
        tags = {'aeroway': True}
        airports = ox.features_from_place(county_name, tags)
        airports = airports[airports.apply(lambda row: 'geometry' in row and row.geometry.geom_type == 'Point', axis=1)]
        airports = airports[airports['aeroway'] == 'aerodrome']
        airports = airports.dropna(subset=['name'])
        airports = airports.sample(min(1, len(airports))) # to make it visible, show 2 random airport
        for idx, row in airports.iterrows():
            if 'geometry' in row and row.geometry.geom_type == 'Point':
                location = (row.geometry.y, row.geometry.x)
                name = row.get('name', 'Unnamed Airport')
                folium.Marker(location=location, popup=name, tooltip=name, icon=folium.Icon(color='darkblue', icon='plane', prefix='fa')).add_to(map_obj)
            
    if show_schools:
        # Let's mark the schools in the county
        tags = {'amenity': 'school'}
        schools = ox.features_from_place(county_name, tags);
        schools = schools[schools.apply(lambda row: 'geometry' in row and row.geometry.geom_type == 'Point', axis=1)]
        schools = schools.dropna(subset=['name'])
        schools = schools.sample(min(3, len(schools))) # to make it visible, show 5 random schools at most
        for idx, row in schools.iterrows():
            if 'geometry' in row and row.geometry.geom_type == 'Point':
                location = (row.geometry.y, row.geometry.x)
                name = row.get('name', 'Unnamed School')
                folium.Marker(location=location, popup=name, tooltip=name, icon=folium.Icon(color='red', icon='plus')).add_to(map_obj)
            
    if show_religious_centers:
        # Let's mark the religious centers in the county
        tags = {'amenity': 'place_of_worship'}
        religious_centers = ox.features_from_place(county_name, tags)
        religious_centers = religious_centers[religious_centers.apply(lambda row: 'geometry' in row and row.geometry.geom_type == 'Point', axis=1)]
        religious_centers = religious_centers.dropna(subset=['name'])
        for idx, row in religious_centers.iterrows():
            if 'geometry' in row and row.geometry.geom_type == 'Point':
                location = (row.geometry.y, row.geometry.x)
                name = row.get('name', 'Unnamed Religious Center')
                folium.Marker(location=location, popup=name, tooltip=name, icon=folium.Icon(color='darkblue', icon='gift')).add_to(map_obj)
                

    return map_obj


def display_map(gdf, show_airports=True, show_schools=True, show_religious_centers=False):
    
    # create a folium map centered on South Carolina
    m = folium.Map(location=[33.8413, -79.0643], zoom_start=7)
    # Add the state boundary with green highlight
    folium.GeoJson(
        gdf,
        style_function=lambda feature: {
            'fillColor': 'gray',
            'color': 'gray',
            'weight': 2,
            'fillOpacity': 0.3
        }
    ).add_to(m)
    
    # let's mark the locations of Fort Bragg and North Auxiliary Airfield on the map
    markers = [
        {"location": [35.139167, -78.999167 ], "name": "Fort Liberty", "color": "green", "icon":"home"}, # ISB
        {"location": [33.616944, -81.083056], "name": "North Auxiliary Airfield", "color": "green", "icon":"home"}, #ISB
        {"location": [33.6504, -84.3737], "name": "FEMA Warehouse", "color": "purple", "icon":"star"}, #FEMA warehouse
        {"location": [34.404063669058125, -81.09810984525734], "name": "SCEMD Warehouse", "color": "green", "icon":"home"}, #SCEMD warehouse
    ] 

    # Add each marker to the map
    for marker in markers:
        folium.Marker(
            location=marker["location"],
            popup=marker["name"],
            tooltip=marker["name"],
            icon=folium.Icon(color=marker["color"], icon=marker["icon"])
        ).add_to(m)
        
    
    for county in counties:
        add_county_to_map(county, m, show_airports, show_schools, show_religious_centers)
        
    return m

m = display_map(gdf, show_airports=True, show_schools=True, show_religious_centers=False)
display(m);

