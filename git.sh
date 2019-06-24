if [ "$2" == "" ];
then
	echo -e "sh $0 <git commit -m > <git add >\nexample: sh $0 'update readme' 'README.md */*' "
	exit
fi
commit=$1
add=$2
uname=`uname|grep MINGW|wc -l`
if [ "$uname" -eq 1 ];
then
	find $add  -type f -exec grep -Iq . {} \; -and -print|xargs -L 1 -I {} dos2unix {}
fi
set -vex
git add $add
git commit -m "$commit"
git push -u origin master

echo done
exit


#create new branch doc(default is in master branch)
git checkout -b doc # will cread new branch doc and switch to doc branch
# do some thing in branch doc then push to remote branch doc by 
git push origin doc
# merge branch doc to master and swith to master branch by 
git checkout master -m 


## for first use git
#git init
#git add *
#git commit -m 'update'
#git remote add origin https://github.com/orangeSi/ClustersPloter.git
#git push -u origin master
#git pull origin master

# config password for git
#git config credential.helper store, store mean forver remember password. git config credential.helper 'cache –timeout=3600' for temperate remember password
# git push orgin master 
# ......
# after this will remember the password

