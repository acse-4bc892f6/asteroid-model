import folium
import webbrowser
from folium import CustomIcon


def plot_circle(lat, lon, radius, map=None, **kwargs):
    """
    Plot a circle on a map (creating a new folium map instance if necessary).

    Parameters
    ----------

    lat: float
        latitude of circle to plot (degrees)
    lon: float
        longitude of circle to plot (degrees)
    radius: float
        radius of circle to plot  (m)
    map: folium.Map
        existing map object

    Returns
    -------

    Folium map object

    Examples
    --------

    >>> import folium
    >>> armageddon.plot_circle(52.79, -2.95, 1e3, map=None)
    """
    map = folium.Map(location=[lat, lon], zoom_start=14)
    map.add_child(folium.LatLngPopup())

    if not map:
        map = folium.Map(location=[lat, lon], control_scale=True)

    folium.Circle([lat, lon], radius, fill=True, fillOpacity=0.6,
                  **kwargs).add_to(map)

    map.save('mapping.html')
    webbrowser.open('mapping.html')

    return map


plot_circle(52.79, -2.95, 1e3, map=None)


def plot_circle_circle(initial_lat, initial_lon, lat, lon, radius, map=None):
    """
    Plot a circle on a map (creating a new folium map instance if necessary).
    Parameters
    ----------
    lat: float
        latitude of circle to plot (degrees)
    lon: float
        longitude of circle to plot (degrees)
    radius: float
        radius of circle to plot (m)
    map: folium.Map
        existing map object
    Returns
    -------
    Folium map object
    Examples
    --------
    >>> import folium
    >>> armageddon.plot_circle(52.79, -2.95, 1e3, map=None)
    """

    map = folium.Map(location=[lat, lon], zoom_start=14)
    map.add_child(folium.LatLngPopup())

    if not map:
        map = folium.Map(location=[lat, lon], control_scale=True)

    initial_location = [initial_lat, initial_lon]
    loction = [lat, lon]

    folium.Marker(location=[lat, lon],
                  popup=loction,
                  icon=folium.Icon(color='red', icon='info-sign')).add_to(map)
    folium.Marker(location=[initial_lat, initial_lon],
                  popup=initial_location,
                  icon=folium.Icon(color='green',
                                   icon='info-sign')).add_to(map)

    url = '{}'.format
    icon_image = url('armageddon/planet.png')
    icon = CustomIcon(icon_image,
                      icon_size=(50, 50),
                      icon_anchor=(22, 22),
                      popup_anchor=(-3, -76))
    folium.Marker(location=[initial_lat, initial_lon], icon=icon,
                  popup='Star').add_to(map)

    folium.Circle(location=[lat, lon],
                  radius=radius[0],
                  color="green",
                  popup="LEVEL1",
                  fill=True).add_to(map)
    folium.Circle(location=[lat, lon],
                  radius=radius[1],
                  color="blue",
                  popup="LEVEL2",
                  fill=True).add_to(map)
    folium.Circle(location=[lat, lon],
                  radius=radius[2],
                  color="yellow",
                  popup="LEVEL3",
                  fill=True).add_to(map)
    folium.Circle(location=[lat, lon],
                  radius=radius[3],
                  color="crimson",
                  popup="LEVEL4",
                  fill=True).add_to(map)

    coordinates = [(initial_lat, initial_lon), (lat, lon)]

    folium.PolyLine(locations=coordinates, weight=2, color='black').add_to(map)

    map.save('mapping.html')
    webbrowser.open('mapping.html')

    return map


def plot_uncertain(uncertain_lat, uncertain_lon, radiu, map=None):
    num = 0
    size = len(uncertain_lat)
    if size == 1:
        num = size
    if size % 2 == 0:
        num = size / 2
    if size % 2 == 0:
        num = (size - 1) / 2

    map = folium.Map(
        location=[uncertain_lat[int(num)], uncertain_lon[int(num)]],
        zoom_start=14)
    map.add_child(folium.LatLngPopup())

    if not map:
        map = folium.Map(
            location=[uncertain_lat[int(size)], uncertain_lon[int(size)]],
            control_scale=True)

    for i in range(len(uncertain_lat)):
        folium.Circle(location=[uncertain_lat[i], uncertain_lon[i]],
                      radius=radiu[i],
                      color="blue",
                      popup="uncertain",
                      fill=True).add_to(map)

    map.save('mapping.html')
    webbrowser.open('mapping.html')
    return map
