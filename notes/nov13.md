

## Analysis of Nov 9 data of Lu-175 hfs scan

```bash
Number of counts recorded: 400957
4712 wavelength recordings with no meaning
Garbage wavelength values: [-2500000.0, -3333333.333]
Duration of operation: 312.5643457
```

In that time, we should expect 312,564 complete buncher cycles. This means that on average, a buncher cycle involves ~1 count. This makes sense, as we had a count rate of about ~1000 hits per second, at a bunching frequency of 1 kHz.


```
Typical wavenumber jitter in a single scan step: 0.0077 cm-1
Typical wavenumber offset from requested: -0.98 cm-1, with scatter 0.02 cm-1
```
 that scatter likely has a large component due to wavelength stability itself.