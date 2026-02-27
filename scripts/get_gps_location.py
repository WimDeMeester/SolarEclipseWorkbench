#!/usr/bin/env python3
"""
Get GPS coordinates for Solar Eclipse Workbench from your smartphone.

Two modes
---------
1. --web  (recommended — no app needed)
   Starts a tiny HTTPS server on this machine.  Open the URL on your phone's
   browser, tap "Get My Location", and the coordinates are sent back here.

   Phone and laptop must be on the same network (WiFi, or use the phone as a
   hotspot — see docs/GPS_PHONE.md for step-by-step instructions).

     python scripts/get_gps_location.py --web

2. --smartphone HOST
   Reads a live NMEA stream from a GPS app running on your phone.

   Recommended Android apps (Google Play, free):
     • "BlueNMEA"  — TCP server, default port 4352
   Recommended iPhone apps (App Store):
     • "GPS2IP Lite" — TCP, default port 11123

     python scripts/get_gps_location.py --smartphone 192.168.1.42
     python scripts/get_gps_location.py --smartphone 192.168.1.42 --smartphone-port 11123

See also: docs/GPS_PHONE.md
"""

import argparse
import socket
import sys
import time


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

FIX_QUALITY = {
    0: "No fix",
    1: "GPS fix",
    2: "DGPS fix (WAAS)",
    3: "PPS fix",
    4: "RTK",
    5: "Float RTK",
    6: "Estimated / dead-reckoning",
}


def format_dms(decimal_degrees: float, is_lat: bool) -> str:
    d = abs(decimal_degrees)
    deg = int(d)
    minutes = (d - deg) * 60
    secs = (minutes - int(minutes)) * 60
    direction = (
        ("N" if decimal_degrees >= 0 else "S")
        if is_lat
        else ("E" if decimal_degrees >= 0 else "W")
    )
    return f"{deg}°{int(minutes):02d}'{secs:05.2f}\"{direction}"


def print_fix(
    lat: float,
    lon: float,
    alt: float | None,
    quality: int = 1,
    num_sats: int = 0,
    hdop: float | None = None,
    speed_knots: float = 0.0,
    timestamp: str = "",
    datestamp: str = "",
) -> None:
    alt_str = f"{alt:.1f} m" if alt is not None else "n/a"
    hdop_str = f"{hdop:.1f}" if hdop is not None else "n/a"
    quality_str = FIX_QUALITY.get(quality, "Unknown")
    print(
        f"\n  Lat: {lat:+.6f}°  ({format_dms(lat, is_lat=True)})\n"
        f"  Lon: {lon:+.6f}°  ({format_dms(lon, is_lat=False)})\n"
        f"  Alt: {alt_str}   Satellites: {num_sats}   HDOP: {hdop_str}\n"
        f"  Fix: {quality_str}   Speed: {speed_knots:.1f} kn\n"
        f"  UTC: {datestamp} {timestamp}\n"
        f"  {'─' * 55}"
    )


def print_copy_paste(lat: float, lon: float, alt: float | None) -> None:
    alt_str = f"{alt:.1f}" if alt is not None else "0.0"
    print("\nReady to paste into Solar Eclipse Workbench:")
    print(f"  Latitude  : {lat:.6f}")
    print(f"  Longitude : {lon:.6f}")
    print(f"  Altitude  : {alt_str} m")


# ---------------------------------------------------------------------------
# Mode 1: browser-based web server (recommended)
# ---------------------------------------------------------------------------

def run_web(port: int = 8765, fix_timeout: float = 300.0) -> None:
    """Start an HTTPS server; open the URL in your phone's browser to submit GPS."""
    try:
        from solareclipseworkbench.phone_gps import WebGpsServer
    except ImportError:
        print("ERROR: solareclipseworkbench package not found.")
        print("  Run from the project root after: pip install -e .")
        sys.exit(1)

    server = WebGpsServer(port=port)

    def on_started(lan_url: str, local_url: str) -> None:
        print(f"\nGPS capture server running ({server.proto.upper()}).\n")
        print("From your PHONE (same WiFi or phone hotspot):")
        print(f"    {lan_url}")
        if server.proto == "https":
            print()
            print("  Certificate warning: tap 'Advanced' → 'Proceed' (Chrome)")
            print("                    or 'Show Details' → 'Visit this website' (Safari)")
        print()
        print("From this LAPTOP (to test):")
        print(f"    {local_url}")
        print()
        print("The page has both automatic GPS and a manual coordinate entry form.")
        print("Server will print a confirmation line when coordinates are received.")
        print("Press Ctrl-C to quit.\n")
        print("  No WiFi at your eclipse site? Use your phone as a hotspot.")
        print("  See docs/GPS_PHONE.md for step-by-step instructions.\n")

    server.start()
    on_started(server.lan_url, server.local_url)

    try:
        fix = server.wait_for_fix(timeout=fix_timeout)
    except KeyboardInterrupt:
        print("\nCancelled by user.")
        server.stop()
        return

    server.stop()

    if fix is None:
        print(f"\nTimeout: no coordinates received within {fix_timeout:.0f}s.")
        return

    print_fix(
        lat=fix["lat"],
        lon=fix["lon"],
        alt=fix.get("alt"),
        num_sats=0,
    )
    print_copy_paste(fix["lat"], fix["lon"], fix.get("alt"))


# ---------------------------------------------------------------------------
# Mode 2: NMEA TCP stream from a smartphone GPS app
# ---------------------------------------------------------------------------

def run_smartphone(host: str, port: int = 4352, fix_timeout: float = 60.0) -> None:
    """Read NMEA sentences streamed from a smartphone GPS app over TCP/WiFi.

    Both phone and laptop must be on the same network.

    Recommended Android apps (Google Play, free):
      • "BlueNMEA"  TCP server, default port 4352

    Recommended iPhone apps (App Store):
      • "GPS2IP Lite"  free, TCP, default port 11123
    """
    try:
        import pynmea2
    except ImportError:
        print("ERROR: pynmea2 is not installed.  Run: pip install pynmea2")
        sys.exit(1)

    print(f"\nConnecting to smartphone GPS at {host}:{port} …")
    print("Make sure your phone and this machine are on the same network.")
    print("Waiting for GPS fix (press Ctrl-C to quit)\n")

    last: dict = {}
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            sock.settimeout(10)
            try:
                sock.connect((host, port))
            except (ConnectionRefusedError, OSError) as exc:
                print(f"\nERROR: Could not connect to {host}:{port} — {exc}")
                print("  • Is the GPS app running and in Server/TCP mode?")
                print("  • Are both devices on the same network?")
                print(f"  • Is the port correct? (default for BlueNMEA: 4352)")
                sys.exit(1)

            print("Connected. Receiving GPS data …\n")
            sock.settimeout(5)
            buf = b""
            deadline = time.monotonic() + fix_timeout
            got_fix = False

            while True:
                try:
                    chunk = sock.recv(4096)
                except socket.timeout:
                    if not got_fix and time.monotonic() > deadline:
                        print(f"\nTimeout: no fix within {fix_timeout:.0f}s.")
                        print("  Make sure the phone has a GPS signal (go outdoors).")
                        return
                    continue

                if not chunk:
                    print("\nConnection closed by phone.")
                    break

                buf += chunk
                while b"\n" in buf:
                    line_bytes, buf = buf.split(b"\n", 1)
                    try:
                        line = line_bytes.decode("ascii", errors="replace").strip()
                    except Exception:
                        continue
                    if not line.startswith("$"):
                        continue
                    try:
                        msg = pynmea2.parse(line)
                    except pynmea2.ParseError:
                        continue

                    if msg.sentence_type == "GGA":
                        quality = int(msg.gps_qual) if msg.gps_qual else 0
                        if quality == 0:
                            elapsed = time.monotonic() - (deadline - fix_timeout)
                            print(f"\r  Searching for satellites … ({elapsed:.0f}s)  ",
                                  end="", flush=True)
                            continue
                        got_fix = True
                        last.update(
                            lat=msg.latitude,
                            lon=msg.longitude,
                            alt=float(msg.altitude) if msg.altitude else None,
                            quality=quality,
                            num_sats=int(msg.num_sats) if msg.num_sats else 0,
                            hdop=float(msg.horizontal_dil) if msg.horizontal_dil else None,
                            speed_knots=last.get("speed_knots", 0.0),
                            timestamp=last.get("timestamp", ""),
                            datestamp=last.get("datestamp", ""),
                        )
                        print("\r", end="")
                        print_fix(**last)

                    elif msg.sentence_type == "RMC" and msg.status == "A":
                        got_fix = True
                        last.update(
                            speed_knots=float(msg.spd_over_grnd) if msg.spd_over_grnd else 0.0,
                            timestamp=str(msg.timestamp) if msg.timestamp else "",
                            datestamp=str(msg.datestamp) if msg.datestamp else "",
                        )

    except KeyboardInterrupt:
        print("\n\nStopped by user.")
        if last:
            print("\nLast recorded fix:")
            print_fix(**last)
            print_copy_paste(last["lat"], last["lon"], last.get("alt"))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Get GPS coordinates for Solar Eclipse Workbench from your smartphone.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--web",
        action="store_true",
        default=False,
        help=(
            "Start a browser-based GPS capture server (recommended, no app needed). "
            "Open the URL shown on your phone's browser and tap 'Get My Location'. "
            "Works over WiFi or phone hotspot."
        ),
    )
    mode.add_argument(
        "--smartphone",
        metavar="HOST",
        default=None,
        help=(
            "Stream NMEA from a smartphone GPS app over TCP "
            "(e.g. --smartphone 192.168.1.42). "
            "Android: 'BlueNMEA' (port 4352). "
            "iPhone: 'GPS2IP Lite' (port 11123)."
        ),
    )
    parser.add_argument(
        "--web-port", type=int, default=8765,
        help="TCP port for --web mode (default: 8765)",
    )
    parser.add_argument(
        "--smartphone-port", type=int, default=4352,
        help="TCP port for --smartphone mode (default: 4352)",
    )
    parser.add_argument(
        "--timeout", type=float, default=300.0,
        help="Seconds to wait for coordinates in --web mode (default: 300)",
    )
    args = parser.parse_args()

    if args.web:
        run_web(port=args.web_port, fix_timeout=args.timeout)
    elif args.smartphone:
        run_smartphone(args.smartphone, port=args.smartphone_port, fix_timeout=args.timeout)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
