def get_git_commit_sha() -> str:
    """
    Retrieve the git commit SHA for the current repository.

    This function uses the `gitpython` library to find the current repository and get the commit SHA of the latest commit.

    Returns:
    - str
        The SHA of the latest git commit.

    Raises:
    - git.InvalidGitRepositoryError
        If the current directory is not part of a git repository.
    - git.NoSuchPathError
        If the path to the git repository is not found.
    """
    import git

    repo = git.Repo(search_parent_directories=True)
    hexsha: str = repo.head.object.hexsha

    return hexsha
