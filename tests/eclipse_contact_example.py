from solareclipseworkbench.solar_eclipse import get_local_circumstances

# Example location: latitude, longitude, and height (meters)
# -3.9852, 41.6669, 828
latitude = 41.6669   # degrees North
longitude = -3.9852  # degrees East (negative for West)
height = 828         # meters above sea level

# latitude = 39.65645555555555555555555555555556
# longitude = -1.466675
# height = 887

result = get_local_circumstances(latitude, longitude, height, "2026-08-12")

def ut_to_hms(ut):
    hours = int(ut)
    minutes = int((ut - hours) * 60)
    seconds = ((ut - hours) * 60 - minutes) * 60
    return hours, minutes, seconds

print("Eclipse Contact Times and Circumstances:")
print(f"Julian Date: {result['jd']}")

max_h, max_m, max_s = ut_to_hms(result['UTMaximum'])
fc_h, fc_m, fc_s = ut_to_hms(result['UTFirstContact'])
sc_h, sc_m, sc_s = ut_to_hms(result['UTSecondContact'])
tc_h, tc_m, tc_s = ut_to_hms(result['UTThirdContact'])
lc_h, lc_m, lc_s = ut_to_hms(result['UTLastContact'])

print(f"Maximum Eclipse UT: {max_h:02}:{max_m:02}:{max_s:05.2f}")
print(f"First Contact UT: {fc_h:02}:{fc_m:02}:{fc_s:05.2f}")
print(f"Second Contact UT: {sc_h:02}:{sc_m:02}:{sc_s:05.2f}")
print(f"Third Contact UT: {tc_h:02}:{tc_m:02}:{tc_s:05.2f}")
print(f"Last Contact UT: {lc_h:02}:{lc_m:02}:{lc_s:05.2f}")

print(f"Eclipse Magnitude: {result['mag']:.6f}")
print(f"Altitude of Sun at Maximum: {result['h']:.2f} degrees")
print(f"Minimum separation (m): {result['m']:.6f}")

duration = (result['UTThirdContact'] - result['UTSecondContact']) * 3600  # seconds
# Convert the duration to minutes and seconds
duration_minutes = int(duration // 60)
duration_seconds = duration % 60
print(f"Totality Duration: {duration_minutes} minutes and {duration_seconds:.2f} seconds")

# Get the next solar eclipses
from solareclipseworkbench.solar_eclipse import get_solar_eclipses

# Get current date
from datetime import datetime, timedelta
current_date = datetime.now()
current_date = current_date.strftime("%Y-%m-%d")  # Format as YYYY-MM-DD
eclipses = get_solar_eclipses(20, current_date)

print("\nUpcoming Solar Eclipses:")
for eclipse in eclipses:
    # Format the date as YYYY-MM-DD with a zero before the month
    print(f"Date: {eclipse['date']}, Type: {eclipse['type']}, Magnitude: {float(eclipse['magnitude']):.6f}, Central duration: {eclipse['duration']}, Saros: {eclipse['saros']}")

