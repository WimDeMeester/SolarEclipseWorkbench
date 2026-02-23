## TODO

### GPS

- [ ] Get GPS coordinates from device
- https://github.com/quentinsf/pygarmin
  
### GUI

### Scripts

### Camera

- [ ] Mirror up - TEST???
- [ ] All commands from Solar Eclipse Maestro: http://xjubier.free.fr/en/site_pages/solar_eclipses/Solar_Eclipse_Maestro_Help/pgs2/btoc6.html

### Problems with gphoto2

- [ ] Camera access only works when executing the command using `sudo`

On macOS, `ptpcamerad` (the system PTP camera daemon) grabs the USB device as soon as gphoto2 releases it, causing `[-53] Could not claim the USB device` errors. There are two permanent fixes:

#### Option A — Image Capture app (recommended, no terminal required, survives macOS updates)

1. Connect the camera
2. Open **Image Capture** (`/Applications/Image Capture.app`)
3. Select the camera in the sidebar
4. Set the bottom-left dropdown **"Connecting this camera opens:"** to **No application**

This tells macOS not to hand the camera to any PTP handler. Persists across reboots.

#### Option B — Disable ptpcamerad permanently via launchctl (macOS 11+)

```bash
sudo launchctl disable system/com.apple.ptpcamerad
sudo launchctl kill TERM system/com.apple.ptpcamerad
```

The `disable` command writes a persistent flag that survives reboots. Undo with:

```bash
sudo launchctl enable system/com.apple.ptpcamerad
```

#### Workaround (temporary, must be repeated after every reboot)

Kill ptpcamerad once per session:
```bash
sudo killall ptpcamerad
```
It will not respawn until the next reboot.