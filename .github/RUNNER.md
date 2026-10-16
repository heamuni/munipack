# CI runner on varda

`ci.yml` runs on a self-hosted GitHub Actions runner on **varda** (Debian 13), set up on 2026-10-14. It is registered to `heamuni/munipack` only, with the name `varda` and the labels `varda` and `docker`; workflows pick it with `runs-on: [self-hosted, varda]`.

| | |
|---|---|
| User | `ghrunner` (system user, no password, in the `docker` group) |
| Runner | `/home/ghrunner/actions-runner` (v2.338.0, updates itself) |
| Job folder | `/home/ghrunner/_work` |
| Service | `actions.runner.heamuni-munipack.varda` (systemd, starts on boot) |
| Docker | `docker.io` and `docker-buildx` from Debian's packages (26.1.5), data in `/home/docker` because `/` is small; build cache capped at 100 GB in `/etc/docker/daemon.json` |

The first build on varda takes longest, since it builds wxWidgets; later builds use the cache. A full run (build, upstream `make check`, `tests/`) took 5 min 20 s on 2026-10-14.

## Checking it

```bash
systemctl status actions.runner.heamuni-munipack.varda      # on varda
```
GitHub: Settings → Actions → Runners should list `varda` as Idle or Active.

## Re-registering

If the runner is removed on GitHub or has to be moved, get a token from Settings → Actions → Runners → New self-hosted runner (it is valid for one hour), then on varda:

```bash
cd /home/ghrunner/actions-runner && sudo ./svc.sh stop && sudo ./svc.sh uninstall
sudo -iu ghrunner bash -c 'cd ~/actions-runner && ./config.sh remove --token TOKEN'
sudo -iu ghrunner bash -c 'cd ~/actions-runner && ./config.sh --unattended --url https://github.com/heamuni/munipack --token TOKEN --name varda --labels varda,docker --work /home/ghrunner/_work'
cd /home/ghrunner/actions-runner && sudo ./svc.sh install ghrunner && sudo ./svc.sh start
```

`svc.sh` writes a systemd unit, so it needs `sudo`; `config.sh` must run as `ghrunner`.
