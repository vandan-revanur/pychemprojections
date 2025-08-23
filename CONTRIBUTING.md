# Contributing

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

We welcome contributions, in the following ways:
1. Issues
2. Enhancement requests (also tracked as GitHub issues.)
3. Pull requests

If you have a question, want to report a bug, or have an enhancement suggestions then please submit an issue.

### Planned Features

We plan to have the following features in the future.

* Sawhorse projection
* Projections/conformations such as Chair/Boat/Twisted boat conformations
* Haworth projections

If you are interested in building any of these join this [discord channel](https://discord.gg/aEmf9MqX)

To contribute with code to the project make a pull request with the below process.

### Pull Request Process

1. Fork this repository.
2. Create a branch: `git checkout -b <branch_name>`.
3. Install pre-commit tool and the hooks: ```pip install pre-commit && pre-commit install```
4. Make your changes and commit them.
Use the [conventional-commits format](https://www.conventionalcommits.org/en/v1.0.0/).
You can do this either with commitizen tool or normal git commit.
Interactive version of the commitizen tool can be found [here](https://github.com/commitizen-tools/commitizen).
For the non-interactive CLI version look [here](https://github.com/streamich/git-cz).
4. [Optionally] Update the README.md with details of changes.
5. [Optionally] Increase the version number and generate changelog using commitizen : `cz bump --changelog`. The versioning scheme we use is [SemVer](http://semver.org/).
4. Push to the remote branch: `git push`
5. Create the pull request.
