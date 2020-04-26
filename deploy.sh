#/bin/bash
git update-index --refresh 
git diff-index --quiet HEAD -- || echo "Untracked changes..."; exit 1
echo "# Before you deploy, have you:"
echo "# 1. Checked that the version number has been bumped, and commited?"
echo "# 2. Checked that the version of software in the README is correct, and commited any changes?"
echo "# If so, run bash deploy.sh | bash"
version=$(cat arbow/version.py | cut -f2 -d "=" | sed -E 's/^[[:space:]]+"(.*)"/v\1/g')
echo "git tag -a $version"
echo "rm -rf build/* dist/*"
echo "python3 setup.py sdist bdist_wheel"
echo "twine upload dist/*"
echo "git push --tags"