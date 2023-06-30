### To be done on your local Mac, linux or WSL machine

### Set an alias to log in quickly

`nano ~/.bashrc` then,  
```bash
alias graham="ssh USERNAME@graham.computecanada.ca"
alias cedar="ssh USERNAME@cedar.computecanada.ca"
```

### Skip entering the password everytime (security warning)

```bash
ssh-keygen
# Go with default name and directory. No passphrase. You have to press enter a few times.
ssh-copy-id -i ~/.ssh/id_rsa USERNAME@CLUSTER.computecanada.ca
```
