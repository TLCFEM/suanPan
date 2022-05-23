try
{
    git | Out-Null
}
catch [System.Management.Automation.CommandNotFoundException]
{
    exit
}

Set-Variable -Name "git_rev" -Value (git rev-parse --short=8 HEAD)

Write-Output "constexpr auto SUANPAN_REVISION = ""$git_rev"";" | Out-File -Encoding utf8 Toolbox/revision.h

echo "Revision tag set to $git_rev"
